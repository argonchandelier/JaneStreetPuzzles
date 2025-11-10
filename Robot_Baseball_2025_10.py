"""
Author:    Logan Richan
Completed: 2025-10-5
Made for solving the Jane Street Puzzle, Oct 2025
"""
import random

class RobotBaseball:
    def __init__(self, p):
        self.p = p
        self.p_inv = 1 - p
        
        """ Notation is: grid[ball #][strike #] """
        self.ks = [[0] * 3 for _ in range(4)]
        self.exp_scores = [[0] * 4 for _ in range(5)]
        self.exp_scores[-1] = [1, 1, 1, 1, -1]
        
        #self.ks_v = [[0] * 3 for _ in range(4)]
        #self.exp_scores_v = [[0] * 4 for _ in range(5)]
        #self.exp_scores_v[-1] = [1, 1, 1, 1, -1]
        
        self.chance_to_get_here = [[0] * 4 for _ in range(5)]
        self.chance_to_get_here[0][0] = 1
    
    def main_process(self):
        if self.p == 0:
            return 0
        
        self.fill_grids()
        self.find_chances_to_arrive()
        
        q = self.chance_to_get_here[3][2]
        
        # print(f"3 balls, 2 strike chance: {q}")
        # print("ks_v:", self.ks_v)
        #print("bbb exp_scores_v:", self.exp_scores_v)
        # print("chance_to_get_here:", self.chance_to_get_here)
        
        return q
    
    def get_k(self, a, c, d):
        num = d - c
        den = a - 2 * c + d
        if num == 0 and den == 0:
            return 1
        
        k = num / den
        return k
    
    def get_exp_score(self, a, c, k):
        return (a - c) * k + c
    
    def fill_grids(self):
        for strik in range(2, -1, -1):
            for ball in range(3, -1, -1):
                c = self.exp_scores[ball][strik + 1]
                d = self.exp_scores[ball + 1][strik]
                a = 4 * self.p + (1 - self.p) * c
                k = self.get_k(a, c, d)
                exp_score = self.get_exp_score(a, c, k)
                
                # print(f"{ball}b, {strik}s; a: {a}, c: {c}, d: {d}, k: {k}, exp_score: {exp_score}")
                
                self.ks[ball][strik] = k
                self.exp_scores[ball][strik] = exp_score
    
    def find_chances_to_arrive(self):
        total_run_chance = 0
        locked_sum = 0
        run_chances = [0] * 6
        for bs in range(0, 6):  # Loop each diagonal
            diagonal_run_sum = 0
            for b in range(bs + 1):
                s = bs - b
                if b > 3 or s > 2:
                    if b == 4 or s == 3:
                        locked_sum += self.chance_to_get_here[b][s]
                    continue
                k = self.ks[b][s]
                max_chance = self.chance_to_get_here[b][s]
                
                ''' Fixed Equations, I got these wrong the first time... '''
                s_chance = ((self.p_inv * k**2) + 2*(1-k)*k)  * max_chance
                b_chance = (1-k)**2 * max_chance
                self.chance_to_get_here[b][s + 1] += s_chance
                self.chance_to_get_here[b + 1][s] += b_chance
                r = (self.p * k**2) * max_chance
                
                diagonal_run_sum += s_chance + b_chance
                total_run_chance += r
                run_chances[bs] += r
            # print(f"+++{bs} diag run sum: {diagonal_run_sum}, total run so far: {total_run_chance},"
            #     f"locked sum: {locked_sum}, total sum: {diagonal_run_sum+total_run_chance+locked_sum}")
        # print(f"Total run chance: {total_run_chance} {run_chances}")


#p0 = 9999/10000
#rb = RobotBaseball(p0)
#print(p0, rb.main_process())


def p_to_q_super(iters=16):
    D = 10
    a, b = 0, 10
    best_ns = []
    best_qs = []
    for iter in range(iters):
        ns = []
        qs = []
        for n in range(a, b + 1):
            p = n/D
            rb = RobotBaseball(p)
            q = rb.main_process()
            
            ns.append(n)
            qs.append(q)
        m = max(qs)
        i = qs.index(m)
        
        bn, bq = ns[i], qs[i]
        best_ns.append(bn)
        best_qs.append(bq)
        
        # print(f"\n++++++++ Best: {bn} / {D}, q: {bq}\n\n")
        print(f"++++++++ Best: q: {bq}, p: {bn} / {D}  (iteration: {iter})")
        
        mid = bn * 10
        a, b = mid - 10, mid + 10
        D *= 10


p_to_q_super()

'''
++++++++ Best: q: 0.2939643532109209, p: 2 / 10  (iteration: 0)
++++++++ Best: q: 0.29594350457911994, p: 23 / 100  (iteration: 1)
++++++++ Best: q: 0.2959679914523962, p: 227 / 1000  (iteration: 2)
++++++++ Best: q: 0.2959679914523962, p: 2270 / 10000  (iteration: 3)
++++++++ Best: q: 0.2959679933462427, p: 22697 / 100000  (iteration: 4)
++++++++ Best: q: 0.29596799337412705, p: 226973 / 1000000  (iteration: 5)
++++++++ Best: q: 0.2959679933742692, p: 2269732 / 10000000  (iteration: 6)
++++++++ Best: q: 0.295967993374272, p: 22697324 / 100000000  (iteration: 7)
++++++++ Best: q: 0.29596799337427215, p: 226973234 / 1000000000  (iteration: 8)
++++++++ Best: q: 0.29596799337427215, p: 2269732337 / 10000000000  (iteration: 9)
++++++++ Best: q: 0.29596799337427215, p: 22697323364 / 100000000000  (iteration: 10)
++++++++ Best: q: 0.29596799337427215, p: 226973233632 / 1000000000000  (iteration: 11)
++++++++ Best: q: 0.29596799337427215, p: 2269732336311 / 10000000000000  (iteration: 12)
++++++++ Best: q: 0.29596799337427215, p: 22697323363103 / 100000000000000  (iteration: 13)
++++++++ Best: q: 0.2959679933742722, p: 226973233631022 / 1000000000000000  (iteration: 14)
++++++++ Best: q: 0.2959679933742722, p: 2269732336310220 / 10000000000000000  (iteration: 15)
'''

class Sim:
    def __init__(self, p, n):
        self.p = p
        self.rb = RobotBaseball(p)
        self.rb.fill_grids()
        #self.rb.find_chances_to_arrive()
        
        self.exp_scores = self.rb.exp_scores
        print("ks:", self.rb.ks)
        print("exp scores:", self.exp_scores)
        
        self.n = int(n)
        self.n_inv = 1/float(n)
        
        self.locations_reached_overall = [[0] * 4 for _ in range(5)]

    def at_bat(self):
        ks = self.rb.ks
        b, s = 0, 0
        q_cond = 0
        locs_reached = [[0] * 4 for _ in range(5)]
        while True: # keep pitching
            locs_reached[b][s] += self.n_inv
            if b == 3 and s == 2:
                q_cond = 1
                
            if s == 3:
                return 0, q_cond, b, s, locs_reached
            if b == 4:
                return 1, q_cond, b, s, locs_reached
            
            k = ks[b][s] # Chance to throw in strike zone
            sw = k # Chance to swing
            rand_k = random.random()
            rand_sw = random.random()
            if rand_k < k: # strikezone pitch
                if rand_sw < sw: # swing
                    rand_hit = random.random()
                    if rand_hit < self.p:
                        locs_reached[4][3] += self.n_inv
                        return 4, q_cond, b, s, locs_reached
                    s += 1
                    continue
                # stay
                s += 1
                continue
            # non-strikezone pitch
            if rand_sw < sw: # swing
                s += 1
                continue
            # stay
            b += 1
            
    def n_sims(self):
        q_conds = 0
        n_inv = 1 / self.n
        for i in range(self.n):
            score, q_cond, b, s, locs_reached = self.at_bat()
            #print(f"score: {score}, q_cond: {q_cond}, b: {b}, s: {s}")
            q_conds += q_cond
            
            for j in range(len(locs_reached)):
                for k in range(len(locs_reached[0])):
                    self.locations_reached_overall[j][k] += locs_reached[j][k]
        
        print(f"locs reached overall: {self.locations_reached_overall}")
        q_est = q_conds / self.n
        print(f"For {self.n} simulations with p={self.p}, q was found to be about {q_est}")

#sim0 = Sim(0.0000001, 100000)
#sim0.n_sims()

def several_sims():
    for i in range(1, 100):
        p = i / 100
        sim = Sim(p, 100000)
        sim.n_sims()

#several_sims()

'''
locs reached overall:   [[0.9999999999999062, 0.7556999999999331, 0.5915999999999512, 0.47989999999996347],
                         [0.24429999999998941, 0.35029999999997774, 0.37479999999997504, 0.2950999999999838],
                         [0.05810000000000064, 0.1310000000000019, 0.17929999999999657, 0.13690000000000124],
                         [0.01419999999999997, 0.04560000000000028, ++0.0880000000000015++, 0.0880000000000015],
                         [0, 0, 0, 0.0001]]
calculated: chance_to_get_here: [[1, 0.4023038618470211, 0.13713128425789378, 0.03427596518584841],
                                 [0.5976559037433533, 0.5620966972196342, 0.3467245598819988, 0.11554788942341375],
                                 [0.30068837285988237, 0.515145299543893, 0.5401616345211114, 0.2699998092169863],
                                 [0.10372301554478063, 0.3097890266880378, ++0.579768979589419++, 0.5794792110070572],
                                 [4.1497499493189634e-05, 0.00012390320920633874, 0.00023181486588941835, 0]]

locs reached overall: [[0.9999999999999062, 0.7618999999999324, 0.5959999999999507, 0.48349999999996307], [0.2380999999999901, 0.3425999999999786, 0.37389999999997514, 0.2924999999999841], [0.06140000000000073, 0.1292000000000021, 0.17769999999999675, 0.13310000000000166], [0.013399999999999974, 0.0463000000000003, 0.09090000000000158, 0.09090000000000158], [0, 0, 0, 0]]

'''
