class B_spline:
    def __init__(self, cp, deg=3):
        self.degree = deg
        self.control_points = cp
        n = self.control_points.shape[1]
        self.knot = [0 for i in range(0,deg)] + np.linspace(0,1,n-deg+1).tolist() + [1 for i in range(0,deg)]
        

    def find_knot(self, u):
        l = len(self.knot)
        if(u == 1):
            return l-self.degree-2
        for a in range(l):
            if(self.knot[a] > u):
                return a-1
        return l


    def evaluate(self, u):
        deg = self.degree
        t0 = self.find_knot(u)
        knot = self.knot[t0-deg:t0+deg+1]
        basis = [0 for t in range(0,deg+1)]
        basis[deg] = 1
        
        for i in range(1,deg+1):
            if(knot[deg-i+1]!=knot[deg+1]):
                basis[deg-i] = basis[deg-i+1]*(knot[deg+1]-u)/(knot[deg+1]-knot[deg-i+1])
            for j in range(deg-i+1,deg):
                if(knot[j]!= knot[j+i]):
                    basis[j] = basis[j]*(u-knot[j])/(knot[j+i]-knot[j])
                if(knot[j+1]!= knot[j+i+1]):
                    basis[j] += basis[j+1]*(knot[j+i+1]-u)/(knot[j+i+1]-knot[j+1])
            if(knot[deg]!=knot[deg+i]):
                basis[deg] *= (u-knot[deg])/(knot[deg+i]-knot[deg])

        val = self.control_points[:,t0-deg:t0+1]
        ans = val.dot(basis)
        return ans


    def plot_data(self,n=100):
        cnt = 0
        data = np.zeros((n+1,3))
        for i in np.linspace(0,1,n+1):
            data[cnt] = self.evaluate(i)
            cnt = cnt+1
        
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot3D(data[:,0],data[:,1], data[:,2])
        # ax.plot3D(data[:,0],data[:,1], data[:,2],'ro')
        ax.plot3D(self.control_points[0,:],self.control_points[1,:],self.control_points[2,:],'ro')

        plt.show()
