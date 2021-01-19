class B_spline_surface:
    def __init__(self, knot_file1, knot_file2, cp_file, deg_file1, deg_file2):
        # self.knot1 = np.genfromtxt(knot_file1, delimiter=' ')
        # self.knot2 = np.genfromtxt(knot_file2, delimiter=' ')
        # self.degree1 = int(np.genfromtxt(deg_file1, delimiter=' '))
        # self.degree2 = int(np.genfromtxt(deg_file2, delimiter=' '))
        # self.control_points = np.genfromtxt(cp_file, delimiter=' ')
        self.knot1 = np.load(knot_file1)
        self.knot2 = np.load(knot_file2)
        self.degree1 = int(np.load(deg_file1))
        self.degree2 = int(np.load(deg_file2))
        self.control_points = np.load(cp_file)
        self.der_requ = 2
        self.der_reqv = 2
        self.reshape(self.knot1.shape[0]+self.degree1+1, self.knot1.shape[0]+self.degree1+1)
        d = self.degree1
        self.knot1 = [0 for i in range(0,d+1)] + self.knot1.tolist() + [1 for i in range(0,d+1)]
        d = self.degree2
        self.knot2 = [0 for i in range(0,d+1)] + self.knot2.tolist() + [1 for i in range(0,d+1)]


    def reshape(self, n, m):
        x = self.control_points
        t0 = np.zeros((n,m))
        t1 = np.zeros((n,m))
        t2 = np.zeros((n,m))
        for i in range(n):
            for j in range(m):
                t0[i,j] = x[i*m+j,0]
                t1[i,j] = x[i*m+j,1]
                t2[i,j] = x[i*m+j,2]
        self.control_points = np.zeros((3,n,m))
        self.control_points[0] = t0
        self.control_points[1] = t1
        self.control_points[2] = t2
        # [t0,t1,t2]
        return 


    def find_knot(self, dir, u):
        if(dir == 'x'):
            knot = self.knot1
            deg = self.degree1
        else:
            knot = self.knot2
            deg = self.degree2
        l = len(knot)
        if(u == 1):
            return l-deg-2
        for a in range(l):
            if(knot[a] > u):
                return a-1
        return l


    def evaluate_basis(self, dir, u):
        if(dir == 'x'):
            deg = self.degree1
            t0 = self.find_knot(dir, u)
            knot = self.knot1[t0-deg:t0+deg+1]
            basis = np.array([0 for t in range(0,deg+1)],dtype=np.float64)
            basis[deg] = 1
        elif(dir == 'y'):
            deg = self.degree2
            t0 = self.find_knot(dir, u)
            knot = self.knot2[t0-deg:t0+deg+1]
            basis = np.array([0.0 for t in range(0,deg+1)],dtype=np.float64)
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
        basis = basis[:,None]
        return [t0, basis]

    def evaluate_basis_derivatives(self, dir, u):
        # direction selection
        if(dir == 'x'):
            deg = self.degree1
            t0 = self.find_knot(dir, u)
            knot = self.knot1[t0-deg:t0+deg+1]
            der_req = self.der_requ;
            # basis = np.array([0 for t in range(0,deg+1)],dtype=np.float64)
            basis = np.zeros((deg+1,deg+1), dtype=np.float64)
            basis[0][deg] = 1
        elif(dir == 'y'):
            deg = self.degree2
            t0 = self.find_knot(dir, u)
            knot = self.knot2[t0-deg:t0+deg+1]
            der_req = self.der_reqv;
            # basis = np.array([0.0 for t in range(0,deg+1)],dtype=np.float64)
            basis = np.zeros((deg+1,deg+1), dtype=np.float64)
            basis[0][deg] = 1
        
        # evaluation of basis upto order degree
        for i in range(1,deg+1):
            if(knot[deg-i+1]!=knot[deg+1]):
                basis[i][deg-i] = basis[i-1][deg-i+1]*(knot[deg+1]-u)/(knot[deg+1]-knot[deg-i+1])
            for j in range(deg-i+1,deg):
                if(knot[j]!= knot[j+i]):
                    basis[i][j] = basis[i-1][j]*(u-knot[j])/(knot[j+i]-knot[j])
                if(knot[j+1]!= knot[j+i+1]):
                    basis[i][j] += basis[i-1][j+1]*(knot[j+i+1]-u)/(knot[j+i+1]-knot[j+1])
            if(knot[deg]!=knot[deg+i]):
                basis[i][deg] = basis[i-1][deg]*(u-knot[deg])/(knot[deg+i]-knot[deg])

        # derivative calculation
        # derivatives = np.zeros((deg+1, der_req+1),dtype=np.float64)
        derivatives = np.zeros((der_req+1, deg+1),dtype=np.float64)

        for i in range(0,deg+1):
            p = deg
            factor = p
            a = np.zeros(der_req+1, dtype=np.float64)
            a[0] = 1
            derivatives[0][i] = basis[deg][i]
            for k in range(1, der_req+1):
                if(i+1 <= deg):
                    deno1 = knot[i+deg+1]-knot[i+k]
                    if(deno1 != 0):
                        a[k] = -a[k-1]/deno1

                for j in range(k-1,0,-1):
                    deno2 = knot[deg+j-k+1+i]-knot[j+i]
                    if(deno2 != 0):
                        a[j] = (a[j]-a[j-1])/deno2
                deno2 = knot[deg-k+1+i]-knot[0+i]
                if(deno2 != 0):
                    a[0] = a[0]/deno2
                rend = min(deg+1, i+k+1)
                # derivatives[i][k] = a[0:rend-i].dot(basis[deg-k][i:rend])
                derivatives[k,i] = factor*a[0:rend-i].dot(basis[deg-k][i:rend])
                p -= 1
                factor *= p

        return [t0,derivatives]

    def plot_data(self,n=100, m=100):
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        ax = plt.axes(projection='3d')  
        cp = self.control_points; d1 = self.degree1; d2 = self.degree2
        new_data = np.zeros((3,self.der_requ+1, self.der_reqv+1,n+1,m+1))
        a = 0; b = 0
        for i in np.linspace(0,1,n+1):
            b = 0
            for j in np.linspace(0,1,n+1):
                # if(a == 60 and b == 80):
                #     pdb.set_trace()
                [x0,basisX] = self.evaluate_basis_derivatives('x',i)
                [y0,basisY] = self.evaluate_basis_derivatives('y',j)
                new_data[:,:,:,a,b] = np.matmul(basisX,np.matmul(cp[:,x0-d1:x0+1,y0-d2:y0+1],basisY.T))
                # new_data[:,:,:,a,b] = np.matmul(basisX.T,np.matmul(cp[:,x0-d1:x0+1,y0-d2:y0+1],basisY))
                b = b+1
            a = a+1
        ax.plot_surface(new_data[0,0,0,:,:],new_data[1,0,0,:,:],new_data[2,0,0,:,:])
        px = 60; py = 80;
        val = new_data[:,0,0,px,py]
        nu = new_data[:,1,0,px,py]; nv = new_data[:,0,1,px,py];
        n = np.cross(nu,nv)
        n /= np.linalg.norm(n)
        nu /= np.linalg.norm(nu)
        nv /= np.linalg.norm(nv)
        line = np.zeros((3,2))
        line[:,0] = val; line[:,1] = val+n;
        ax.plot(line[0,:], line[1,:], line[2,:])
        u,v = np.meshgrid(np.linspace(-1,1,101),np.linspace(-1,1,101))
        ax.plot_surface(u+val[0],v+val[1],(-u*n[0]-v*n[1])/n[2]+val[2])
        
        plt.show()
