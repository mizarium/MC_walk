import random
import numpy as np
import matplotlib.pyplot as plt
from tqdm.notebook import tqdm, trange

class data():
    
    def __init__(self, lengths, starts, ends, target):
        self.lengths = lengths
        self.starts = starts
        self.ends = ends
        self.target = target
        self.len = len(lengths)


class MC_walk():
    
    def __init__ (self, n_reps=10, stop_val=10, max_iter=10):
        self.data = ''
        self.n_reps = n_reps
        self.stop = stop_val
        self.max_iter = max_iter
        self.Poff = np.linspace (0, 1, 6)
        self.Pend = np.linspace (0, 1, 6)
        self.SE = np.zeros((6, 6))
        self.Poff_log = []
        self.Pend_log = []
        self.SE_log = []
        self.poff_best = []
        self.pend_best = []
        self.pwin_predicted = ''
        
        
    def walk (self, l, start, end, p_off, p_end):
        
        # this functin generates one walk along the substrate with 
        # length l, starting point 'start', end poind 'end'
        
        stop_criterion = False
        position = start
        
        while position != end:
                
            x = random.random()
            y = random.random()
                    
            if y <= p_off:
                break                   
            
            elif position == 1:
                if y <= p_off + p_end:
                    break
                else:
                    position = position + 1
                    
            elif position == l:
                if y <= p_off + p_end:
                    break
                else:
                    position = position - 1
            
            elif x > 0.5:
                position = position + 1
            else:
                position = position - 1
        
        if position == end:
            return 1
        else:
            return 0
        
        
    def walk10000 (self, l, start, end, p_off, p_end, nwalks=100, nlaunches=100):
        
        # This functin generates 1000 walks along the substrate with 
        # length l, 
        # starting point 'start'
        # end poind 'end'
        # 
        # and returns # of successfull walks over each 100 
    
        result = []
        for _ in range(nlaunches):
            success = 0
            fail = 0
            
            for __ in range(nwalks):
                walk_result = self.walk(l, start, end, p_off, p_end)
                success += walk_result
                fail += 1 - walk_result

            result.append(success)
            assert success + fail == 100
        
        return np.array(result)
        
        
    def calculate_cum_error (self, array, target):
        
        array = array/100
        return np.sum((array - target)**2)
        
        
    def fit(self, data):
        
        self.data = data
        for _ in range (self.n_reps):
            
            self.Poff = np.linspace (0, 1, 6)
            self.Pend = np.linspace (0, 1, 6)
            
            SE_iter_log = []
            Poff_iter_log = [self.Poff]
            Pend_iter_log = [self.Pend]
            
            for __ in range (self.max_iter):
                
                # SE - matrix, containing error values
                
                SE = np.zeros((6, 6))
                
                for i in trange(len(self.Poff)):
                    for j in range(len(self.Pend)):
                        SE_next = 0
                        
                        for k in range (self.data.len):
                            
                            # simulating walks in both directions
                            
                            w1 = self.walk10000(data.lengths[k], self.data.starts[k], self.data.ends[k], self.Poff[i], self.Pend[j])
                            w2 = self.walk10000(data.lengths[k], self.data.ends[k], self.data.starts[k], self.Poff[i], self.Pend[j])
                            w_err1 = self.calculate_cum_error (w1, self.data.target[k])
                            w_err2 = self.calculate_cum_error (w2, self.data.target[k])
                            
                            SE_next += (w_err1 + w_err2)
                            
                        SE[i, j] = SE_next
                
                # fit function tries to converge to the minimum error point
                
                best = np.where(SE == np.min(SE))
                criterion = np.std(SE) < np.min(SE)/10 
                        
                poff_best = self.Poff[best[0]]
                pend_best = self.Pend[best[1]]
                
                #update logs
                SE_iter_log.append(SE)
                Poff_iter_log.append(self.Poff)
                Pend_iter_log.append(self.Pend)
                        
                if criterion: 
                    print ('Yahhhoo')
                    break
                    
                # if criterion is not reached, the boundries for next modeling step 
                # are defined by following rule:
                
                if best[0] == 0:
                    next_poff_left = self.Poff[0]
                    next_poff_right = self.Poff[1]
                elif best[0] == 5:
                    next_poff_left = self.Poff[-2]
                    next_poff_right = self.Poff[-1]
                else:
                    next_poff_left = self.Poff[int(best[0]) - 1]
                    next_poff_right = self.Poff[int(best[0]) + 1]

                if best[1] == 0:
                    next_pend_left = self.Pend[0]
                    next_pend_right = self.Pend[1]        
                elif best[1] == 5:
                    next_pend_left = self.Pend[-2]
                    next_pend_right = self.Pend[-1]
                else:
                    next_pend_left = self.Pend[int(best[1]) - 1]
                    next_pend_right = self.Pend[int(best[1]) + 1]
       
                self.Poff =  np.linspace (next_poff_left, next_poff_right, 6)
                self.Pend =  np.linspace (next_pend_left, next_pend_right, 6)       
       
                print ('next Poff = ', self.Poff)
                print ('next Pend = ', self.Pend)
            
            
            # update logs
            self.SE_log.append(SE_iter_log)
            self.Poff_log.append(Poff_iter_log)
            self.Pend_log.append(Pend_iter_log)
            self.poff_best.append(np.mean(self.Poff))
            self.pend_best.append(np.mean(self.Pend))
        
        self.evaluate()
            
            
    def evaluate(self):
        
        # when modeling is finished, MC_walk estimates probability of a successfull 
        # random walk (the enzyme reaches 'end' point) given mean values of p_off and p_end
        # obtained by n_reps iterations
        
        poff_final = np.mean(self.poff_best)
        pend_final = np.mean(self.pend_best)
        prediction = []
        
        for i in range(self.data.len):
            pwin = self.walk10000(self.data.lengths[i], self.data.starts[i], self.data.ends[i], poff_final, pend_final) #cahnge here for back walk
            prediction.append(np.mean(pwin)/100)
        
        self.pwin_predicted = prediction
        print (prediction)

    
    def out_round(self, num):
        if num > 1:
            return round(num, 1)
        elif num > 0.01:
            return round(num, 3)
        elif num == 0:
            return 0
        else:
            return '%.1e' % num
        
        
    def plot_iteration(self, num, verbose=True):
        
        # this function helps to visualise shape of the objective function by iteration
        # if verbose == True, actual values of estimated error are shown
        
        LEN = len(self.SE_log[num])
        if LEN <= 8:
            nrows = 2
        else:
            nrows = 3

        fig, ax = plt.subplots(nrows, 4, figsize=(18 , nrows*4))
        fig.add_subplot(111, frameon=False)
        
        for i in range(LEN):
            x = i//4
            y = i%4
            ax[x, y].imshow(self.SE_log[num][i])
            ax[x, y].set_xticks(np.arange(6))
            ax[x, y].set_yticks(np.arange(6))
            ax[x, y].set_xticklabels([self.out_round(j) for j in self.Pend_log[num][i+1]], fontsize=10, position=(0, 0.03))
            ax[x, y].set_yticklabels([self.out_round(j) for j in self.Poff_log[num][i+1]], fontsize=10, position=(0.03, 0))
            plt.setp(ax[x, y].get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
   
            if verbose:
                for a in range(6):
                    for b in range(6):
                        text = ax[x, y].text(b, a, self.out_round(self.SE_log[num][i][a, b]),
                                   ha="center", va="center", color="lightgrey")
    
        for i in range(LEN, nrows*4, 1):
            ax[i//4, i%4].set_xticks([])
            ax[i//4, i%4].set_yticks([])
        
        plt.xticks([])
        plt.yticks([])
        plt.xlabel('p end', fontsize = 20, labelpad=30)
        plt.ylabel('p off', fontsize = 20, labelpad=30)
        plt.title ('Error by iteration',  fontsize = 25, pad=25)


if __name__ == "__main__":
         
    Lengths = [41, 61, 81, 101]
    Starts = [8, 8, 8, 8]
    Ends = [29, 49, 69, 89]
    T = [0.3925, 0.2643, 0.1933, 0.1199]

    DATA = data(Lengths, Starts, Ends, T)
    my_walk = MC_walk(n_reps=10, max_iter=10)
    my_walk.fit(DATA)
    poff_final = np.mean(my_walk.poff_best)
    pend_final = np.mean(my_walk.pend_best)
    print ("")
    print ("##############################################")
    print ("")
    print ("Modeling results: ")
    print ("p_off = ", poff_final)
    print ("p_end = ", pend_final)
    print ("Predicted Pwin for given data")
    print (my_walk.prediction)
    #my_walk.plot_iteration(0)
