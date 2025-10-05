"""
Author:    Logan Richan
Completed: 2025-10-5
Made for solving the Jane Street Puzzle, Oct 2025
"""
class RobotBaseball:
    def __init__(self, p):
        self.p = p
        self.p_inv = 1 - p
        
        """ Notation is: grid[ball #][strike #] """
        self.ks = [[0] * 3 for _ in range(4)]
        self.exp_scores = [[0] * 4 for _ in range(5)]
        self.exp_scores[-1] = [1, 1, 1, 1, -1]
        
        self.ks_v = [[0] * 3 for _ in range(4)]
        self.exp_scores_v = [[0] * 4 for _ in range(5)]
        self.exp_scores_v[-1] = [1, 1, 1, 1, -1]
        
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
        # print("exp_scores_v:", self.exp_scores_v)
        # print("chance_to_get_here:", self.chance_to_get_here)
        
        return q
    
    def get_k(self, a, c, d):
        num = d - c
        den = a - 2 * c + d
        if (num) == 0 and (den) == 0:
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
                self.ks_v[ball][strik] = k
                self.exp_scores_v[ball][strik] = exp_score
    
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
                # print(b, s)
                k = self.ks_v[b][s]
                max_chance = self.chance_to_get_here[b][s]
                
                # print(f"pinv: {self.p_inv}")
                s_chance = (self.p_inv * k) * max_chance
                b_chance = (1 - k) * max_chance
                self.chance_to_get_here[b][s + 1] += s_chance
                self.chance_to_get_here[b + 1][s] += b_chance
                r = (self.p * k) * max_chance
                # print(s_chance, b_chance, r)
                
                diagonal_run_sum += s_chance + b_chance
                total_run_chance += r
                run_chances[bs] += r
            # print(f"+++{bs} diag run sum: {diagonal_run_sum}, total run so far: {total_run_chance},"
            #     f"locked sum: {locked_sum}, total sum: {diagonal_run_sum+total_run_chance+locked_sum}")
        # print(f"Total run chance: {total_run_chance} {run_chances}")


p0 = 9999/10000
rb = RobotBaseball(p0)
print(p0, rb.main_process())


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
++++++++ Best: q: 0.4072225223084541, p: 1 / 10  (iteration: 0)
++++++++ Best: q: 0.5651223931250657, p: 1 / 100  (iteration: 1)
++++++++ Best: q: 0.5784789260601996, p: 1 / 1000  (iteration: 2)
++++++++ Best: q: 0.579768979589419,  p: 1 / 10000  (iteration: 3)
++++++++ Best: q: 0.5798975042384327, p: 1 / 100000  (iteration: 4)
++++++++ Best: q: 0.5799103518712396, p: 1 / 1000000  (iteration: 5)
++++++++ Best: q: 0.5799116365861743, p: 1 / 10000000  (iteration: 6)
++++++++ Best: q: 0.5799117650571841, p: 1 / 100000000  (iteration: 7)
++++++++ Best: q: 0.5799117779042803, p: 1 / 1000000000  (iteration: 8)
++++++++ Best: q: 0.5799117791889898, p: 1 / 10000000000  (iteration: 9)
++++++++ Best: q: 0.5799117793174609, p: 1 / 100000000000  (iteration: 10)
++++++++ Best: q: 0.5799117793303079, p: 1 / 1000000000000  (iteration: 11)
++++++++ Best: q: 0.5799117793315927, p: 1 / 10000000000000  (iteration: 12)
++++++++ Best: q: 0.5799117793317212, p: 1 / 100000000000000  (iteration: 13)
++++++++ Best: q: 0.5799117793317341, p: 1 / 1000000000000000  (iteration: 14)
++++++++ Best: q: 0.5799117793317351, p: 1 / 10000000000000000  (iteration: 15)

Ans: q: 0.5799117793
'''



