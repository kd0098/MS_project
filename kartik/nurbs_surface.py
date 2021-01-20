class nurbs_surface(B_spline_surface):
    def __init__(self, cp, weight, n, m, deru=2, derv=2, deg1=3, deg2=3):
        B_spline_surface.__init__(self, cp, n, m, deru, derv, deg1, deg2)
        self.weight = weight.reshape(self.control_points[0].shape)
        self.control_points *= self.weight
        self.control_points = np.append(self.control_points, [self.weight], axis = 0)
    

    def plot_data(self,n=100, m=100):
        cp = self.control_points; d1 = self.degree1; d2 = self.degree2
        new_data = np.zeros((4,self.der_requ+1, self.der_reqv+1,n+1,m+1))
        a = 0; b = 0
        for i in np.linspace(0,1,n+1):
            b = 0
            for j in np.linspace(0,1,n+1):
                [x0,basisX] = self.evaluate_basis_derivatives('x',i)
                [y0,basisY] = self.evaluate_basis_derivatives('y',j)
                new_data[:,:,:,a,b] = np.matmul(basisX,np.matmul(cp[:,x0-d1:x0+1,y0-d2:y0+1],basisY.T))
                b = b+1
            a = a+1

        # store combinations for fast operations
        ncr = np.zeros((self.der_requ+1, self.der_reqv+1))
        for i in range(1,self.der_requ+1):
            ncr[i,1] = 1
            for j in range(2,i+1):
                ncr[i,j] = ncr[i-1,j-1]+ncr[i-1,j]

        # evaluation of derivatives of nurbs including function value
        for i in range(self.der_requ+1):
            for j in range(self.der_reqv+1):
                for t in range(1,i):
                    new_data[0:3,i,j,:,:] -= ncr[i,t]*new_data[0:3,i-t,j,:,:]*new_data[3,t,0,:,:]
                for s in range(1,j):
                    new_data[0:3,i,j,:,:] -= ncr[j,s]*new_data[0:3,i,j-s,:,:]*new_data[3,0,s,:,:]
                for t in range(1,i):
                    for s in range(1,j):
                        new_data[0:3,i,j,:,:] -= ncr[i,t]*ncr[j,s]*new_data[0:3,i-t,j-s,:,:]*new_data[3,t,s,:,:]
                new_data[0:3,i,j,:,:] /= new_data[3,0,0,:,:]

        # surface
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        ax = plt.axes(projection='3d')
        ax.plot_surface(new_data[0,0,0,:,:],new_data[1,0,0,:,:],new_data[2,0,0,:,:])
        # return new_data
        
        plt.show()
