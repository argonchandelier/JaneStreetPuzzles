"""
Author:    Logan Richan
Completed: 2025-09-24
"""

'''
Dividing the Problem:
xMake all combos of Ls by corner in>out or out>in; only 4 of each total: l,r,t,b and one of t,b and one of l,r each
xCheck all combos of #s in each L for each configuration above based on given numbers (including outside box); some will have 0 possibilities
xOf valid of last step, Utilize given pentomino shapes to go through all possibilities (location based, from outside of box) that will:
    -Have a value sum in each pentomino that adds to a multiple of 5
    -Not go over the limit of # of #s in each L
    -(Not cross other placed pentomino)
    -(Not violate the 1st-pentomino-called-from-outside rule)
    -...Check every spot from origin outward in each rotation/reflection
    -..."Grow" the pentomino square by square
x(Eliminate violators of 2x2-regions-have-at-least-1-empty-square rule)
    x? Check board and locking empty squares (as empty) that have any diagonal + other 2 adjacent-to-both squares (4 checks each empty sq))
xFor any boards valid so far, find any pentominos that fill in the rest (choose 3 of 6 possibilities from (FLPTWY)):
    -(Remember invalid squares given by outside grid letters)
    -Scout every hole to check for all possibilities to even place any of the 6 choices (check mult-of-5 & L-#-of-#s rules each space)
    -If at least 3 different ones can be placed, check all combinations of valid placements (checking 2x2 rule)
    -(...Potentially a lot of checks, but a lot has been eliminated up to now, so hopefully just a (relatively) few # of possibilities,
         especially as there can only be 1 valid board after this)
    -? 2x2 rule final check last
xFinally, there should only be 1 solution, so use a final answer calculator

Also consider:
-Seeing number of possible boards after certain steps
-Have a board visualizer, for doing manual checks; perhaps several ways of checking too

-Last FLPTWY pentominos: Search for 2x2 empty areas and have attempts to place each of these there...
'''
import copy

class Grid:
    def __init__(self):
        self.grid = [[0] * 9 for _ in range(9)] # '0' is empty or unfilled
        self.regionSize = [[0] * 9 for _ in range(9)] # Based on size, not # inside
        self.regionConvertType = [0] * 9 # index 0 is size 1, etc.; number stored is the actual number appearing in the region; 0 means undetermined
        self.blanks = [[0] * 9 for _ in range(9)] # 1s are guaranteed blank; 2 is filled; 0 is undetermined; 3 is pre-filled; 4: used pre-filled (in pentomino)
        
        self.row_of_3 = -1
        self.row_of_7 = -1
        self.col_of_6 = -1
        self.col_of_2 = -1

        self.grid[0][4] = 5
        self.grid[1][3] = 4
        self.grid[4][4] = 1
        self.grid[7][5] = 8
        self.grid[8][4] = 9

        self.blanks[0][4] = 3
        self.blanks[1][3] = 3
        self.blanks[4][4] = 3
        self.blanks[7][5] = 3
        self.blanks[8][4] = 3

        self.regionSize[4][4] = 1
        self.regionConvertType[0] = 1

    def SetSpace(self, row, col, num, fourOK = False, specialSpace = False):
        blankIndex = self.GetBlanksSpace(row, col)
        if num == 0:
            self.grid[row][col] = 0
            self.blanks[row][col] = 0
            return True
        if blankIndex == 1 or (blankIndex == 4 and not fourOK):
            return False

        self.grid[row][col] = num
        if blankIndex == 3: # prefilled number has special condition (only og 5 #s)
            self.blanks[row][col] = 4
        elif blankIndex == 4:
            self.blanks[row][col] = 3
        elif specialSpace:
            self.blanks[row][col] = 3
        else:
            self.blanks[row][col] = 2

        return True

    def ResetSpace(self, row, col):
        self.grid[row][col] = 0

        if self.blanks[row][col] == 4:
            self.blanks[row][col] = 3
        else:
            self.blanks[row][col] = 0

    def SetSpaceToRegionTypeTry(self, row, col):
        if row < 0 or row > 8 or col < 0 or col > 8:
            return -1
        blankType = self.GetBlanksSpace(row, col)
        if blankType > 0 and blankType != 3:
            return -1

        regType = self.GetRegionType(row, col)
        self.SetSpace(row, col, regType)
        return regType

    def GetSpace(self, row, col):
        return self.grid[row][col]

    def GetBlanksSpace(self, row, col):
        return self.blanks[row][col]
    
    def SetBlanksSpace(self, row, col, index):
        self.blanks[row][col] = index

    def ApplyRegionSize(self, row, col, size):
        self.regionSize[row][col] = size

        # Only applied early on anyway, so might as well have here
        if self.grid[row][col] > 0:
            self.SetRegionType2(size, self.grid[row][col])

    def SetRegionType(self, row, col, typeNum):
        size = self.regionSize[row][col]
        self.regionConvertType[size - 1] = typeNum

    def SetRegionType2(self, size, typeNum):
        self.regionConvertType[size - 1] = typeNum

    def RegionSizeToType(self, size):
        return self.regionConvertType[size-1]

    def GetRegionType(self, row, col):
        size = self.regionSize[row][col]
        return self.RegionSizeToType(size)

    def SetRegionToTypeIfValidOuterCheck(self, row, col, typeNum):
        regType = self.GetRegionType(row, col)
        if regType == 0:
            self.SetRegionType(row, col, typeNum)
            self.SetSpace(row, col, typeNum, specialSpace=True)
            return True
        return False

    def ResetRegionOuterCheck(self, row, col):
        self.SetRegionType(row, col, 0)
        self.SetSpace(row, col, 0)

    def DisplayRegionSize(self):
        s = "Regions (size #s):\n"
        for r in range(9):
            s += str(self.regionSize[r]) + "\n"
        print(s)

    def DisplayGrid(self):
        s = "Grid:\n"
        for r in range(9):
            s += str(self.grid[r]) + "\n"
        print(s)

    def DisplayBlanks(self):
        s = "Guaranteed grid blanks:\n"
        for r in range(9):
            s += str(self.blanks[r]) + "\n"
        print(s)

    def DisplayAll(self, gridNo=""):
        print(f"\n----------GRID {gridNo}\n")
        self.DisplayGrid()
        self.DisplayRegionSize()
        self.DisplayBlanks()
        print(str(self.regionConvertType))
        print("\n----------END\n")

    def Rule2x2Blank(self):
        for r in range(8):
            for c in range(8):
                summ = 0
                space = []
                for r2 in [r, r+1]:
                    for c2 in [c, c+1]:
                        if self.GetSpace(r2, c2):
                            summ += 1
                        else:
                            space = [r2, c2]
                if summ == 4:
                    return False
                elif summ == 3:
                    self.blanks[space[0]][space[1]] = 1
        return True

    def RuleNNumsInRegionTypeNum(self):
        counts = [0]*9
        for r in range(9):
            for c in range(9):
                n = self.grid[r][c]
                if n > 0:
                    size = self.regionSize[r][c]
                    counts[size-1] += 1

        for i, N in enumerate(self.regionConvertType):
            count = counts[i]
            if count > N:
                return False
            # If count == N, then all other in region would be blank=1
        return True
    
    def RuleRecheckOuterNums(self):
        for c in range(0, 9):
            s = self.GetSpace(3, c)
            if s > 0:
                if s != 6:
                    return False
                break
        for c in range(8, -1 , -1):
            s = self.GetSpace(5, c)
            if s > 0:
                if s != 2:
                    return False
                break
        for r in range(0, 9):
            s = self.GetSpace(r, 6)
            if s > 0:
                if s != 7:
                    return False
                break
        for r in range(8, -1, -1):
            s = self.GetSpace(r, 2)
            if s > 0:
                if s != 3:
                    return False
                break
        return True
    
    def RuleFinalNoLooseNums(self):
        # All nums are part of pentomino, so check for mult of 5
        count = 0
        for r in range(9):
            for c in range(9):
                s = self.GetSpace(r, c)
                if s > 0:
                    count += 1
        
        if count % 5 == 0:
            return True
        return False
    
    def RuleFinalConnectedRegion(self):
        filled_spaces = []
        chains = []
        for r in range(9):
            for c in range(9):
                regionType1 = self.GetBlanksSpace(r, c)
                if regionType1 < 2:
                    continue
                
                origin_space = [r, c]
                if origin_space in filled_spaces:
                    continue
                filled_spaces.append([r, c])
                
                # Start a chain
                chain = [origin_space]
                latest_spaces = [origin_space]
                while True:
                    next_latest_spaces = []
                    for space in latest_spaces:
                        for add_space in [[0, 1], [0, -1], [1, 0], [-1, 0]]:
                            r2, c2 = space[0] + add_space[0], space[1] + add_space[1]
                            pot_new_space = [r2, c2]
                            if r2 < 0 or r2 > 8 or c2 < 0 or c2 > 8:
                                continue
                            if pot_new_space in filled_spaces:
                                continue
                            
                            regionType2 = self.GetBlanksSpace(r2, c2)
                            if regionType2 < 2:
                                continue
                            
                            filled_spaces.append(pot_new_space)
                            next_latest_spaces.append(pot_new_space)
                            chain.append(pot_new_space)
                    latest_spaces = next_latest_spaces
                    if len(next_latest_spaces) == 0:
                        break
                
                chains.append(chain)
        
        n_chains = len(chains)
        
        if n_chains == 1:
            return True
        return False
            
    def FindFinalAns(self):
        empty_spaces = []
        '''
        groups = []
        for r in range(9):
            for c in range(9):
                if self.grid[r][c] > 0:
                    empty_spaces.append([r, c])
                    groups.append([[r,c]])
        
        groupNum = 0
        while True:
            cdsInGroup = 0
            while True:
                origGroups = copy.deepcopy(groups) #groups[groupNum][cdsInGroup]
                '''
        chains = []
        for r in range(9):
            for c in range(9):
                #if self.grid[r][c] == 0:
                #    continue
                regionType1 = self.GetBlanksSpace(r, c)
                if regionType1 > 1:
                    continue
                
                origin_space = [r,c]
                if origin_space in empty_spaces:
                    continue
                empty_spaces.append([r, c])
                
                # Start a chain
                chain = [origin_space]
                latest_spaces = [origin_space]
                #print(f"origin space: {origin_space}")
                while True:
                    #print("latest spaces:", latest_spaces)
                    next_latest_spaces = []
                    for space in latest_spaces:
                        for add_space in [[0, 1], [0, -1], [1, 0], [-1, 0]]:
                            r2, c2 = space[0]+add_space[0], space[1]+add_space[1]
                            pot_new_space = [r2, c2]
                            if r2 < 0 or r2 > 8 or c2 < 0 or c2 > 8:
                                continue
                            if pot_new_space in empty_spaces:
                                continue
                            
                            regionType2 = self.GetBlanksSpace(r2, c2)
                            if regionType2 > 1:
                                continue
                            
                            empty_spaces.append(pot_new_space)
                            #print(empty_spaces)
                            next_latest_spaces.append(pot_new_space)
                            chain.append(pot_new_space)
                    latest_spaces = next_latest_spaces
                    if len(next_latest_spaces) == 0:
                        break
                        
                chains.append(chain)
        
        product = 1
        for chain in chains:
            n = len(chain)
            product *= n
        #print(f"n empty spaces: {len(empty_spaces)}, {empty_spaces}")
        
        return product
        



class Hooks:
    def __init__(self):
        self.grids = []
        self.debugCount = 0
        
        # For last 3 pentomino placements
        self.validPentsPlaced = -1
        self.valid2x2s = []
        
        self.currentGrid = Grid()
        

    def MainProcess(self):
        self.FindAllRegionCombos()

        self.UseNumbersOutsideBox()
        print("After Outer Numbers; Len of new self.grids:", len(self.grids))  # Full: 68358
        #self.grids[0].DisplayAll(0) if len(self.grids) > 0 else print("No grids")

        self.PentominoStart()
        print(f"After Outer Pentominos; num new Grids: {len(self.grids)}") # Full: 78300 (old:10150)
        #self.grids[0].DisplayAll(0) if len(self.grids) > 0 else print("No grids")
        #self.grids[100].DisplayAll(100) if len(self.grids) > 100 else print("...")
        
        self.Rules2Negate()
        print(f"After Rules2Negate #1; num new Grids: {len(self.grids)}") # Full: 1003
        #self.grids[0].DisplayAll(0) if len(self.grids) > 0 else print("No grids")
        
        self.FinalPentominos()
        print(f"After Final Pentominos; num new Grids: {len(self.grids)}")  # Full: 3978
        #self.grids[0].DisplayAll(0) if len(self.grids) > 0 else print("No grids")
        
        self.Rules2Negate()
        print(f"After Rules2Negate #2; num new Grids: {len(self.grids)}")  # Full: 8
        #self.grids[0].DisplayAll(0) if len(self.grids) > 0 else print("No grids")
        #for i, grid in enumerate(self.grids):
        #    grid.DisplayAll(i)
        
        self.FinalCheck()
        print(f"After FinalCheck; num new Grids: {len(self.grids)}")  # Full: 2 (They are the same)
        for i, grid in enumerate(self.grids):
            grid.DisplayAll(i)

    def FindAllRegionCombos(self):
        def fillLRegion(lrDir, tbDir, row, col, layer, grid):
            size = 9-layer

            for c in range(col, col + size*lrDir, lrDir):
                grid.ApplyRegionSize(row, c, size)
            for r in range(row, row + size*tbDir, tbDir):
                grid.ApplyRegionSize(r, col, size)

        def deeper(t, b, l, r, layer, grid):
            if layer == 8:
                self.grids.append(grid)
                return

            if t > 0 and l > 0:
                newGrid = copy.deepcopy(grid)
                fillLRegion(1, 1, 4-t, 4-l, layer, newGrid)
                deeper(t-1, b, l-1, r, layer+1, newGrid)
            if t > 0 and r > 0:
                newGrid = copy.deepcopy(grid)
                fillLRegion(-1, 1, 4 - t, 4 + r, layer, newGrid)
                deeper(t - 1, b, l, r - 1, layer + 1, newGrid)
            if b > 0 and l > 0:
                newGrid = copy.deepcopy(grid)
                fillLRegion(1, -1, 4 + b, 4 - l, layer, newGrid)
                deeper(t, b-1, l-1, r, layer+1, newGrid)
            if b > 0 and r > 0:
                newGrid = copy.deepcopy(grid)
                fillLRegion(-1, -1, 4 + b, 4 + r, layer, newGrid)
                deeper(t, b-1, l, r-1, layer+1, newGrid)

        deeper(4, 4, 4, 4, 0, Grid())
        # Every permutation of 4ts+4ds combined with every permutation of 4ls+4rs: ( 8!/(4!*4!) )**2 = 4900
        print(len(self.grids)) # displayed: 4900; checks out
        #self.grids[1811].DisplayRegionSize()

    def UseNumbersOutsideBox(self):
        newGrids = []
        #for grid in [self.grids[0]]:
        for grid in self.grids:
            def placeOuter6():
                for c1 in range(0, 9):
                    # Valid is if L region is not already for a different #
                    valid = grid.SetRegionToTypeIfValidOuterCheck(3, c1, 6)
                    if c1 > 0:
                        grid.blanks[3][c1-1] = 1
                    if valid:
                        grid.col_of_6 = c1
                        placeOuter2()
                        grid.ResetRegionOuterCheck(3, c1)
                for c in range(0, 9):
                    grid.blanks[3][c] = 0

            def placeOuter2():
                for c2 in range(8, -1, -1):
                    valid = grid.SetRegionToTypeIfValidOuterCheck(5, c2, 2)
                    if c2 < 8:
                        grid.blanks[5][c2+1] = 1
                    if valid:
                        grid.col_of_2 = c2
                        placeOuter3()
                        grid.ResetRegionOuterCheck(5, c2)
                for c in range(8, -1, -1):
                    grid.blanks[5][c] = 0

            def placeOuter3():
                onRow = 8
                for r1 in range(8, -1, -1):
                    if grid.GetSpace(r1, 2) > 0 or grid.GetBlanksSpace(r1, 2) == 1: # 6 or 2 made it first...
                        break
                    onRow = r1
                    valid = grid.SetRegionToTypeIfValidOuterCheck(r1, 2, 3)
                    if r1 < 8:
                        grid.blanks[r1+1][2] = 1
                    if valid:
                        grid.row_of_3 = r1
                        placeOuter7()
                        grid.ResetRegionOuterCheck(r1, 2)
                for r1 in range(8, onRow, -1):
                    grid.blanks[r1][2] = 0

            def placeOuter7():
                onRow = -1
                for r2 in range(0, 9):
                    if grid.GetSpace(r2, 6) > 0 or grid.GetBlanksSpace(r2, 6) == 1: # 6 or 2 made it first...
                        break
                    onRow = r2
                    valid = grid.SetRegionToTypeIfValidOuterCheck(r2, 6, 7)
                    if r2 > 0:
                        grid.blanks[r2-1][6] = 1
                    if valid:
                        grid.row_of_7 = r2
                        newGrid = copy.deepcopy(grid)
                        newGrids.append(newGrid)
                        grid.ResetRegionOuterCheck(r2, 6)
                for r2 in range(0, onRow):
                    grid.blanks[r2][6] = 0

            placeOuter6()
        self.grids = newGrids
        

    def PentominoStart(self):
        newGrids = []
        
        base_region_N = [[0, 0], [-1, 0], [-2, 0], [-2, 1], [-3, 1]]

        for grid in self.grids:
            def fillRegions(act_regions, next_func, R, dir):
                for act_region in act_regions:
                    summ = 0
                    break_i = -1
                    for i, coords in enumerate(act_region):
                        r = coords[0]
                        c = coords[1]
                        n = grid.SetSpaceToRegionTypeTry(r, c)
                        if n == -1:
                            summ = -1
                            break_i = i
                            break
                        summ += n
                    if summ % 5 == 0:
                        ### Turn blank guaranteed here
                        start = 0 if dir == 1 else 8
                        end = 8 if dir == 1 else 0
                        new_blank_cs = []
                        for c in range(start, end, dir):
                            blankIndex = grid.GetBlanksSpace(R, c)
                            if blankIndex == 2 or blankIndex == 4:
                                break
                            if blankIndex == 0:
                                grid.SetBlanksSpace(R, c, 1)
                                new_blank_cs.append(c)
                        ###
                        next_func()
                        ### Turn blank guaranteed back here
                        for new_blank_c in new_blank_cs:
                            grid.SetBlanksSpace(R, new_blank_c, 0)
                        ###
                    for i, coords in enumerate(act_region): # Revert back
                        r, c, n = coords[0], coords[1], coords[2]
                        if i == break_i:
                            break
                        grid.SetSpace(r, c, n, fourOK=True)

            def placeI():
                rel_region0 = [[0,0], [0,1], [0,2], [0,3], [0,4]] # Horiz
                rel_region1 = [[0,0], [1,0], [2,0], [3,0], [4,0]] # Vert
                act_regions = []
                R = 0
                for C in range(0, 3): # Horiz
                    r = 0
                    act_region = []
                    for coords in rel_region0:
                        c = C + coords[1]
                        act_region.append([r, c, grid.GetSpace(r, c)])
                    act_regions.append(act_region)
                for C in range(0, 7): # Vert
                    act_region = []
                    for coords in rel_region1:
                        r = R + coords[0]
                        c = C + coords[1]
                        act_region.append([r, c, grid.GetSpace(r, c)])
                    act_regions.append(act_region)

                fillRegions(act_regions, placeU, R, 1)
                '''
                for act_region in act_regions:
                    sum = 0
                    for coords in act_region:
                        r = coords[0]
                        c = coords[1]
                        n = grid.SetSpaceToRegionType(r, c)
                        sum += n
                    if sum % 5 == 0:
                        continue
                    placeU()
                    for coords in act_region: # Revert back
                        r, c, n = coords[0], coords[1], coords[2]
                        grid.SetSpace(r, c, n)
                        '''

            def placeU():
                rel_region0 = [[0, 0], [1, 0], [1, -1], [1, -2], [0, -2]] # u
                rel_region1 = [[0, 0], [0, -1], [1, -1], [2, -1], [2, 0]] # c
                rel_region2 = [[0, 0], [1, 0], [0, -1], [0, -2], [1, -2]] # n
                rel_region3 = [[0, 0], [0, -1], [1, 0], [2, 0], [2, -1]] # D
                act_regions = []
                R = 0
                ''' Remember to check for invalid squares like with the 7s '''
                endC = 1
                if grid.row_of_7 == R:
                    endC = 5
                for rel_region in [rel_region0, rel_region1, rel_region2, rel_region3]:
                    for C in range(8, endC, -1):
                        act_region = []
                        #print(rel_region)
                        for coords in rel_region:
                            r, c = R+coords[0], C+coords[1]
                            #print(r, c)
                            act_region.append([r, c, grid.GetSpace(r, c)])
                        act_regions.append(act_region)

                fillRegions(act_regions, placeZ, R, -1)

            def placeZ():
                rel_regionS = [[0,0], [0, 1], [-1, 1], [-2, 1], [-2, 2]]
                rel_regionSp = [[-2, 0], [-1, 0], [-1, 1], [-1, 2], [0, 2]]
                rel_regionZ = [[-2, 0], [-2, 1], [-1, 1], [0, 1], [0, 2]]
                rel_regionZp = [[0, 0], [-1, 0], [-1, 1], [-1, 2], [-2, 2]]
                act_regions = []
                R = 8
                endC = 7
                
                for rel_region in [rel_regionS, rel_regionSp, rel_regionZ, rel_regionZp]:
                    if grid.row_of_3 == R:
                        endC = 3
                        if rel_region == rel_regionZ:
                            endC = 2
                        if rel_region == rel_regionSp:
                            endC = 1
                    for C in range(0, endC):
                        act_region = []
                        for coords in rel_region:
                            r, c = R+coords[0], C+coords[1]
                            act_region.append([r, c, grid.GetSpace(r, c)])
                        act_regions.append(act_region)
                
                fillRegions(act_regions, placeV, R, 1)

            def placeV():
                #print("V")
                rel_regionBR = [[0, 0], [0, -1], [0, -2], [-1, 0], [-2, 0]]
                rel_regionBL = [[0,0], [0,-1], [0,-2], [-1,-2], [-2,-2]]
                rel_regionTL = [[-2,-2], [-1,-2], [0,-2], [-2,-1], [-2,0]]
                rel_regionTR = [[0,0], [-1,0], [-2,0], [-2,-1], [-2,-2]]
                act_regions = []
                R = 8
                for rel_region in [rel_regionBR, rel_regionBL, rel_regionTL, rel_regionTR]:
                    for C in range(8, 2, -1):
                        act_region = []
                        for coords in rel_region:
                            r, c = R + coords[0], C + coords[1]
                            act_region.append([r, c, grid.GetSpace(r, c)])
                        act_regions.append(act_region)
                
                fillRegions(act_regions, placeX, R, -1)

            def placeX():
                rel_regionU = [[0, -1], [-1, -1], [-1, 0], [-2, -1], [-1, -2]]
                rel_regionM = [[0, 0], [0, -1], [0, -2], [-1, -1], [1, -1]]
                rel_regionD = [[1, 0], [1, -1], [1, -2], [0, -1], [2, -1]]
                act_regions = []
                R = 3
                endC = 1
                for rel_region in [rel_regionU, rel_regionM, rel_regionD]:
                    if grid.row_of_7 == R:
                        endC = 5
                        if rel_region == rel_regionM:
                            endC = 4
                    for C in range(8, endC, -1):
                        act_region = []
                        for coords in rel_region:
                            r, c = R + coords[0], C + coords[1]
                            act_region.append([r, c, grid.GetSpace(r, c)])
                        act_regions.append(act_region)
                
                fillRegions(act_regions, placeN, R, -1)
        
            def placeN():
                
                R = 5
                base_region = [[0, 0], [-1, 0], [-2, 0], [-2, 1], [-3, 1]]
                '''
                Spread/Grow Method: 1 grow = 2-stick, add 3-stick to each of those + 2 grow (3-stick) with 2-stick to each of those
                2: udr; 3: umdr;; 2u-3,2d-3; 2r-4; 3u-3,3m-4,3d-3; 3r-4 = 24
                '''
                reg_2utl = [[0, 0], [-1, 0], [-1, -1], [-2, -1], [-3, -1]]
                reg_2utr = [[0, 0], [-1, 0], [-1, 1], [-2, 1], [-3, 1]]
                reg_2ubr = [[0, 0], [-1, 0], [0, 1], [1, 1], [2, 1]]
                reg_2dtr = [[0, 0], [1, 0], [0, 1], [-1, 1], [-1, 2]]
                reg_2dbl = [[0, 0], [1, 0], [1, -1], [2, -1], [3, -1]]
                reg_2dbr = [[0, 0], [1, 0], [1, 1], [2, 1], [3, 1]]
                reg_2rtl = [[0, 0], [0, 1], [-1, 0], [-1, -1], [-1, 2]]
                reg_2rtr = [[0, 0], [0, 1], [-1, 1], [-1, 2], [-1, 3]]
                reg_2rbl = [[0, 0], [0, 1], [1, 0], [1, -1], [1, -2]]
                reg_2rbr = [[0, 0], [0, 1], [1, 1], [1, 2], [1, 3]]
                reg_3utl = [[0, 0], [-1, 0], [-2, 0], [-2, -1], [-3, -1]]
                reg_3utr = [[0, 0], [-1, 0], [-2, 0], [-2, 1], [-3, 1]]
                reg_3ubr = [[0, 0], [-1, 0], [-2, 0], [0, 1], [1, 1]]
                reg_3mtl = [[-1, 0], [0, 0], [1, 0], [-1, -1], [-2, -1]]
                reg_3mtr = [[-1, 0], [0, 0], [1, 0], [-1, 1], [-2, 1]]
                reg_3mbl = [[-1, 0], [0, 0], [1, 0], [1, -1], [2, -1]]
                reg_3mbr = [[-1, 0], [0, 0], [1, 0], [1, 1], [2, 1]]
                reg_3dtr = [[0, 0], [1, 0], [2, 0], [0, 1], [-1, 1]]
                reg_3dbl = [[0, 0], [1, 0], [2, 0], [2, -1], [3, -1]]
                reg_3dbr = [[0, 0], [1, 0], [2, 0], [2, 1], [3, 1]]
                reg_3rtl = [[0, 0], [0, 1], [0, 2], [-1, 0], [-1, -1]]
                reg_3rtr = [[0, 0], [0, 1], [0, 2], [-1, 2], [-1, 3]]
                reg_3rbl = [[0, 0], [0, 1], [0, 2], [1, 0], [1, -1]]
                reg_3rbr = [[0, 0], [0, 1], [0, 2], [1, 2], [1, 3]]
                rel_regions = [reg_3dbl, reg_3dbr, reg_3dtr, reg_3mbr, reg_3rtl, reg_3mtl, reg_3rbr, reg_3rbl, reg_2dbl, reg_2dbr, reg_2dtr, reg_2rbl, reg_2rbr,
                               reg_3rtr, reg_2rtr, reg_2rtl, reg_2ubr, reg_3ubr, reg_2utr, reg_3utr, reg_3mbl, reg_3mtr, reg_3utl, reg_2utl]
                act_regs = []
                
                endC = 9
                if grid.row_of_3 == R:
                    endC = 3
                for rel_region in rel_regions:
                    for C in range(0, endC):
                        act_region = []
                        broken = False
                        for coords in rel_region:
                            r, c = R + coords[0], C + coords[1]
                            if r < 0 or r > 8 or c < 0 or c > 8:
                                broken = True
                                break
                            act_region.append([r, c, grid.GetSpace(r, c)])
                        if not broken:
                            act_regs.append(act_region)
                
                fillRegions(act_regs, finalOutside, R, 1)

            def finalOutside():
                newGrid = copy.deepcopy(grid)
                newGrids.append(newGrid)
                #newGrid.DisplayAll(self.debugCount)
                self.debugCount += 1

            placeI()

        self.grids = newGrids
        
    def Rules2Negate(self):
        newGrids = []
        for grid in self.grids:
            valid = True
            valid = grid.Rule2x2Blank()
            if not valid:
                continue
            # 461
            valid = grid.RuleNNumsInRegionTypeNum()
            if not valid:
                continue
            
            valid = grid.RuleRecheckOuterNums()
            
            if valid:
                newGrids.append(grid)
                # 0???
                
        self.grids = newGrids
    
    def FinalPentominos(self):
        newGrids = []
        
        # Left: FLPTWY
        # L:8; P:8; T:8; W:4; Y:8; F:8
        # F note: Have the center 2x2 where one of the squares outside is diagonal to the apex / the 2 outer squares are adjacent
        def checkPenominos(onFunc):
            '''
            -Check if the 4/8 outer spots are 0 or 3 blank
            -Ensure sum is 5
            -Place and move to next pent in list
            '''
            self.validPentsPlaced += 1
            if self.validPentsPlaced == 3:
                newGrid = copy.deepcopy(self.currentGrid)
                newGrids.append(newGrid)
                self.validPentsPlaced -= 1
                return
            
            funcs = [placeL, placeP, placeT, placeW, placeY, placeF]
            for i, func in enumerate(funcs[onFunc : 4 + self.validPentsPlaced]):
                func()
            
            self.validPentsPlaced -= 1
        
        def fillRegion(func_num, rel_regions):
            for sq in self.valid2x2s:
                R, C = sq[0], sq[1]
                for rel_region in rel_regions:
                    act_region = []
                    broken = False
                    for coords in rel_region:
                        r, c = coords[0] + R, coords[1] + C
                        if r < 0 or r > 8 or c < 0 or c > 8:
                            broken = True
                            break
                        act_region.append([r, c, self.currentGrid.GetSpace(r, c)])
                    if broken:
                        continue
                    
                    summ = 0
                    break_i = -1
                    for i, coords in enumerate(act_region):
                        r = coords[0]
                        c = coords[1]
                        n = self.currentGrid.SetSpaceToRegionTypeTry(r, c)
                        if n == -1:
                            summ = -1
                            break_i = i
                            break
                        summ += n
                    if summ % 5 == 0:
                        #next_func()
                        checkPenominos(func_num+1)
                    for i, coords in enumerate(act_region):  # Revert back
                        r, c, n = coords[0], coords[1], coords[2]
                        if i == break_i:
                            break
                        self.currentGrid.SetSpace(r, c, n, fourOK=True)
        
        # These functions go through all valid2x2s (already found) and place all 4/8 on each one
        # If valid found, call checkPenominos(nextIndex)
        def placeL():
            rel0 = [[0,0], [1,0], [0,1], [0,2], [0,3]]
            rel1 = [[0,0], [0,1], [1,0], [2,0], [3,0]]
            rel2 = [[0,-2], [0,-1], [0,0], [0,1], [1,1]]
            rel3 = [[0,0], [0,1], [1,1], [2,1], [3,1]]
            rel4 = [[1,0], [1,1], [0,1], [-1,1], [-2,1]]
            rel5 = [[0,1], [1,1], [1,0], [1,-1], [1,-2]]
            rel6 = [[-2,0], [-1,0], [0,0], [1,0], [1,1]]
            rel7 = [[0,0], [1,0], [1,1], [1,2], [1,3]]
            rel_regions = [rel0, rel1, rel2, rel3, rel4, rel5, rel6, rel7]
            
            fillRegion(0, rel_regions)
            
        def placeP():
            rel0 = [[0,0], [0,1], [1,0], [1,1], [-1, 0]]
            rel1 = [[0,0], [0,1], [1,0], [1,1], [-1, 1]]
            rel2 = [[0,0], [0,1], [1,0], [1,1], [0, -1]]
            rel3 = [[0,0], [0,1], [1,0], [1,1], [1, -1]]
            rel4 = [[0,0], [0,1], [1,0], [1,1], [0, 2]]
            rel5 = [[0,0], [0,1], [1,0], [1,1], [1, 2]]
            rel6 = [[0,0], [0,1], [1,0], [1,1], [2, 0]]
            rel7 = [[0,0], [0,1], [1,0], [1,1], [2, 1]]
            rel_regions = [rel0, rel1, rel2, rel3, rel4, rel5, rel6, rel7]
            
            fillRegion(1, rel_regions)
        
        def placeT():
            rel0 = [[0,0], [0,1], [1,0], [0,-1], [2,0]]
            rel1 = [[0,0], [0,1], [1,0], [-1,0], [0,2]]
            rel2 = [[0,0], [0,1], [1,1], [0,2], [-2,1]]
            rel3 = [[0,0], [0,1], [1,1], [0,-1], [-1,1]]
            rel4 = [[0,0], [1,0], [1,1], [2,0], [1,2]]
            rel5 = [[0,0], [1,0], [1,1], [-1,0], [1,-1]]
            rel6 = [[0,1], [1,0], [1,1], [2,1], [1,-1]]
            rel7 = [[0,1], [1,0], [1,1], [1,2], [-1,1]]
            rel_regions = [rel0, rel1, rel2, rel3, rel4, rel5, rel6, rel7]
            
            fillRegion(2, rel_regions)
        
        def placeW():
            rel0 = [[0,0], [0,1], [1,0], [-1,1], [1,-1]]
            rel1 = [[0,0], [0,1], [1,1], [-1,0], [1,2]]
            rel2 = [[0,0], [1,0], [1,1], [0,-1], [2,1]]
            rel3 = [[0,1], [1,0], [1,1], [0,2], [2,0]]
            rel_regions = [rel0, rel1, rel2, rel3]
            
            fillRegion(3, rel_regions)
        
        def placeY():
            rel0 = [[0,0], [0,1], [1,0], [0,-1], [0,-2]]
            rel1 = [[0,0], [0,1], [1,0], [-1,0], [-2,0]]
            rel2 = [[0,0], [0,1], [1,1], [-1,1], [-2,1]]
            rel3 = [[0,0], [0,1], [1,1], [0,2], [0,3]]
            rel4 = [[0,0], [1,0], [1,1], [2,0], [3,0]]
            rel5 = [[0,0], [1,0], [1,1], [1,-1], [1,-2]]
            rel6 = [[0,1], [1,0], [1,1], [2,1], [3,1]]
            rel7 = [[0,1], [1,0], [1,1], [1,2], [1,3]]
            rel_regions = [rel0, rel1, rel2, rel3, rel4, rel5, rel6, rel7]
            
            fillRegion(4, rel_regions)
        
        def placeF():
            rel0 = [[0,0], [0,1], [1,0], [-1,-1], [-1,0]]
            rel1 = [[0,0], [0,1], [1,0], [-1,-1], [0,-1]]
            rel2 = [[0,0], [0,1], [1,1], [-1,2], [0,2]]
            rel3 = [[0,0], [0,1], [1,1], [-1,2], [-1,1]]
            rel4 = [[0,0], [1,0], [1,1], [2,-1], [0,2]]
            rel5 = [[0,0], [1,0], [1,1], [2,-1], [1,-1]]
            rel6 = [[0,1], [1,0], [1,1], [2,2], [1,2]]
            rel7 = [[0,1], [1,0], [1,1], [2,2], [2,1]]
            rel_regions = [rel0, rel1, rel2, rel3, rel4, rel5, rel6, rel7]
            
            fillRegion(5, rel_regions)
        
        self.valid2x2s = []
        for r in range(8):
            for c in range(8):
                valid_spaces = 0
                for r2 in [r, r + 1]:
                    for c2 in [c, c + 1]:
                        stat = self.currentGrid.GetBlanksSpace(r2, c2)
                        if stat == 0 or stat == 3:
                            valid_spaces += 1
                if valid_spaces == 4:
                    self.valid2x2s.append([r, c])
        
        for gridx in self.grids:
            self.currentGrid = gridx
            checkPenominos(0)
            
        self.grids = newGrids

    def FinalCheck(self):
        newGrids = []
        
        good_count = 0
        for grid in (self.grids):
            valid = grid.RuleFinalNoLooseNums()
            if not valid:
                continue
            valid = grid.RuleFinalConnectedRegion()
            if not valid:
                continue
                
            newGrids.append(grid)
            ans = grid.FindFinalAns()
            print(f"Grid {good_count} Answer: {ans}")
            good_count += 1
        
        self.grids = newGrids


hooks11 = Hooks()
hooks11.MainProcess()

