import networkx 
import matplotlib.pyplot as plt
import numpy as np
import collections
import math
import mpmath



def nba_graph(n, m, kernel=None, initial_graph=None, seed=None):
 k=1 
 j=0
 ds=[0]

 if m < 1 or m >= n:
     raise networkx.NetworkXError(
        f"Barabási–Albert network must have m >= 1 and m < n, m = {m}, n = {n}"
     )

 if kernel is None:
     def kernel(x):
         return x    

 if initial_graph is None:
     G = networkx.star_graph(m)
 else:
        if len(initial_graph) < m or len(initial_graph) > n:
            raise networkx.NetworkXError(
                f"initial graph needs between m={m} and n={n} nodes"
            )
        G = initial_graph.copy() 
 
 ds[0]=G.degree[0]

 while k<=m: 
     ds.append(G.degree[k])  
     k += 1
 
 source = len(G)
      
 while source < n:
        dist = [kernel(d) for d in ds]
        targets = networkx.utils.discrete_sequence(m, distribution=dist, seed=seed)
        G.add_edges_from(zip([source]*m, targets))

        while j<m:
            ds[targets[j]] += 1
            j+=1
        
        ds.append(m)  
        source += 1
        j=0
 return G


o=nba_graph(10000,1,kernel=lambda x: x ** (-0.5))


degree_freq_o = networkx.degree_histogram(o)
degrees_o = range(len(degree_freq_o))

c=1
#c=3

#mu=2*c #gamma=+1
#mu=1.32724 #gamma=+0.5, c=1
#mu=2.34124 #gamma=+0.5, c=3
mu=0.79049 #gamma=-0.5, c=1
#mu=0.440216 #gamma=-0.5, c=3
#mu=0.66815 #gamma=-0.9, c=1
#mu=0.230863 #gamma=-0.9, c=3
#mu=1 #gamma=0


#A= 1.107449907 #gamma=+1, c=1
#A= 2.596384714  #gamma=+1, c=3
#A=25.52874473 #gamma=+0.5, c=1
#A=16.13957783 #gamma=+0.5, c=3
A=21.9765828 #gamma=-0.5, c=1
#A=6198935.48 #gamma=-0.9, c=1


# приближено разпределение за gamma=1 - в края

#1>gamma>0, приближено:
#z = np.linspace(4, max(degrees_o))
#q=A*mu/(c+mu/pow(c,0.5))*pow(z,(pow(mu/c,2)-1)/2)*pow(math.e,-2*mu/c*pow(z,0.5))*10000 #gamma=+0.5

#-1<gamma<0, приближено:
q=[]
z=[]
for l in range(2,max(degrees_o)+10):
    #gamma=-0.5
    qz=10000*A*mu/c/(mu/c+pow(c,(-0.5)))/pow(mu/c,l-c)*pow(math.factorial(c-1),0.5)/pow(math.factorial(l-1),0.5)*pow(math.e,-pow((mu/c),-1)/0.5*pow(l,0.5)+0.5*pow(mu/c,-2)*np.log(l))
    #gamma=-0.9
    #qz= 10000*A*mu/c/(mu/c+pow(c,(-0.9)))/pow(mu/c,l-c)*pow(math.factorial(c-1),0.9)/pow(math.factorial(l-1),0.9)*pow(math.e,-c/mu/0.1*pow(l,0.1))
    
    q.append(qz)
    l=l+1
    z.append(l)



#gamma=+1, #точно:
#z = np.linspace(4, max(degrees_o)) 
#r=2*c*(c+1)/z/(z+1)/(z+2)*10000 

#gamma=+0.5, gamma=-0.5, gamma=-0.9, gamma=0; точно:
#r=mu/c/pow(c,0.5)/(1+mu/c/pow(c,0.5))*10000
r=mu/c/pow(c,-0.5)/(1+mu/c/pow(c,-0.5))*10000
#r=mu/c/pow(c,-0.9)/(1+mu/c/pow(c,-0.9))*10000 
#r=mu/c/pow(c,0)/(1+mu/c/pow(c,0))*10000 
rs=[r]
zs=[c]
k=c+1
for k in range (c+1,56):
    #r=r*pow(k-1,0.5)/pow(k,0.5)/(1+mu/c/pow(k,0.5))    
    r=r*pow(k-1,-0.5)/pow(k,-0.5)/(1+mu/c/pow(k,-0.5))
    #r=r*pow(k-1,-0.9)/pow(k,-0.9)/(1+mu/c/pow(k,-0.9))
    #r=r*pow(k-1,0)/pow(k,0)/(1+mu/c/pow(k,0))
    rs.append(r)
    zs.append(k)
      


fig = plt.figure(figsize=(8, 8))
axgrid = fig.add_gridspec(ncols=2, nrows=1,  width_ratios=[1,1]) #дава възможност в едно изображение да се съдържат мрежата и разпределението ѝ

ax0 = fig.add_subplot(axgrid[0,0])
#networkx.draw(o, with_labels = False, node_size=22 , node_color="darkcyan", width=0.3)

ax2 = fig.add_subplot(axgrid[0,1])
ax2.scatter(degrees_o[1:], degree_freq_o[1:],label='генерирано разпределение', color='darkcyan' ,clip_on=False) #bar(*np.unique(degree_sequence, return_counts=True))
#ax2.plot(z,r, label='точно разпределение', color='r') #gamma=+1
ax2.plot(zs,rs, label='точно разпределение', color='r') #gamma=+0.5, gamma=-0.5, gamma=-0.9

#gamma=+1, приближено:
#z = np.linspace(4, max(degrees_o)) 
#q=mu/(c+(mu/pow(c,1)))*A*pow(z,-3)*10000 

ax2.plot(z, q, label='приближено разпределение', color='black') #всички приближени разпределения

plt.ylim(1,None)

ax2.set_xscale("log")
ax2.set_yscale("log");
ax2.set_xlabel("Степен")
ax2.set_ylabel("# върхове")
ax2.legend()

fig.tight_layout()
plt.show()