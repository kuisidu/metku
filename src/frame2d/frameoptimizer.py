import numpy as np
#from frame2d import Frame2D, LineLoad
from scipy.optimize import brute, basinhopping, dual_annealing
from timeit import default_timer as timer

class FrameOptimizer:
    
    def __init__(self,
                 frame,
                 opt_algorithm='brute',
                 lb=None,
                 ub=None,
                 costf=None,
                 constf=None,
                 maxiter=100):
    
        if not frame.is_generated:
            raise TypeError('Frame must be generated before optimization')
        self.frame = frame
        if len(lb):
            self.lb = lb
        else:
            self.lb = [50] * len(frame.members)           
        if len(ub):
            self.ub = ub
        else:
            self.ub = [1000] * len(frame.members)  
        if len(self.lb) != len(self.ub):
            raise ValueError('Lower and upper boundary must be same length')
        self.opt_algorithm = opt_algorithm
        self.costf = costf
        self.constf = constf
        self.maxiter = maxiter
        
      
    def cost_function(self, X, ratio_limit):
        
        R = self.constraint_function(X, ratio_limit)
        R = np.clip(R, -100, 0)
        return self.frame.weight + sum(R**2 * 1e5)
    
    def constraint_function(self, X, ratio_limit=1):
        
        self.change_values(X)
        
        R = []
        for member in self.frame.members.values():
            R.extend(member.r)
            
        return ratio_limit - np.asarray(R)
     
        
    def change_values(self, X):
        for i, member in enumerate(self.frame.members.values()):
            val = round(X[i])
            member.h = val
            #print(member.profile)
        self.frame.calculate()
        
    def optimize(self, ratio_limit=1, debug=True):
        """ Optimizes frame with given optimization algorithm
        
            Returns:
            --------
            :return frame: Optimized frame
            :rtype frame: Frame2D
        """
        
        start = timer()
        
        args = ratio_limit
        if not self.constf:
            constraint_function = self.constraint_function
        else:
            constraint_function = self.constf
        if not self.costf:
            cost_function = self.cost_function
        else:
            cost_function = self.costf
        # initial guess
        x0=[]
        bounds = []
        for i in range(len(self.lb)):
            bounds.append((self.lb[i], self.ub[i]))

        for i in bounds:
            x0.append(np.random.uniform(i[0],i[1]))
        

           
        if debug:
            print("Initial guess x0: ", x0)
            iprint = 2
        else:
            iprint = 0
       
        if self.opt_algorithm == 'brute':
            res = brute(cost_function,
                        bounds,
                        disp=debug)
            print("Optimal values: ", res)
            
        elif self.opt_algorithm == 'basinhopping':
            kwargs = {'args': args}
            res = basinhopping(cost_function,
                                x0,
                                minimizer_kwargs=kwargs,
                                stepsize=30,
                                niter=self.maxiter,
                                interval=10,
                                disp=debug)
            self.change_values(res.x)
            print("Optimal values: ", res.x)
        
        elif self.opt_algorithm == 'dual_annealing':
            pass
        end = timer()
        print(f'Time elapsed: {end - start :.2f} s')
        return self.frame
    


if __name__ == '__main__':
    
    frame = Frame2D(simple=[2,1,5,10], supports='fixed', num_elements=2)
    frame.add(LineLoad(frame.members[4], [-50, -50], 'y'))
    frame.add(LineLoad(frame.members[5], [-50, -50], 'y'))
    #frame.add(LineLoad(frame.members[8], [-50, -50], 'y'))
    #frame.add(LineLoad(frame.members[9], [-50, -50], 'y'))
    frame.hinge_joints()
    frame.generate()
    optimizer = FrameOptimizer(frame,
                               'basinhopping',
                               maxiter=200)
    opt_frame = optimizer.optimize(debug=True, ratio_limit=1)
    opt_frame.plot()
    
