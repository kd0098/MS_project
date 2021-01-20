class beizer_curve:
    def __init__(self, cp):
        self.control_points = cp

    def evaluate(self,u):
        temp = self.control_points.copy()
        for i in range(len(temp.columns)):
            for j in range(len(temp.columns)-i-1):
                temp.iloc[:,j] = (1-u)*temp.iloc[:,j]+u*temp.iloc[:,j+1]
        temp = temp.iloc[:,0]
        return temp

    def plot_data(self,n=100):
        data = pd.DataFrame()
        a = 0
        for i in np.linspace(0,1,n+1):
            data[a] = self.evaluate(i)
            a = a+1
        
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot3D(data.iloc[0,:],data.iloc[1,:], data.iloc[2,:])
        ax.plot3D(self.control_points.iloc[0,:],self.control_points.iloc[1,:],self.control_points.iloc[2,:],'ro')
        plt.show()
