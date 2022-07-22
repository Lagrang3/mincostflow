import networkx 
import random
import sys
import tempfile
import numpy
import subprocess
import matplotlib.pyplot as plt

def generate(n,m,mcap,mwei,f):
    S=0
    T=1
    g=networkx.gnm_random_graph(n,m,directed=True)
    print(n,m,S,T,file=f)
    for a,b in g.edges():
        print(a,b,random.choice(range(mcap)),random.choice(range(mwei)),file=f)

def plot_data(series):
    fig = plt.figure()
    ax = fig.subplots()
    for x in series:
        ax.loglog(n_array,series[x],'-o',label=x)
    ax.legend()
    ax.set_xlabel('N nodes')
    ax.set_ylabel('Time (micro seconds)')
    fig.savefig('latest.png')

if __name__ == "__main__":
    n_array = [ 2**i for i in range(7,13)]
    m_array = [ int(n*7.5) for n in n_array ]
    cost_arr = [ 200 ]*len(n_array)
    cap_arr = [ 200 ]*len(n_array)
    rep_arr = [20]*len(n_array)
    
    names = ['Ortools','Edmonds-Karp','Primal-dual','Capacity-scaling']
   
    series = dict()
    for x in names:
        series[x] = numpy.zeros(len(n_array))
    
    for i in range(len(n_array)):
        n,m,cap,cost = n_array[i],m_array[i],cap_arr[i],cost_arr[i]
        print("N nodes = ",n)        
        for j in range(rep_arr[i]):
            print("   rep = ",j)
            fname = tempfile.mktemp()
            with open(fname,'w') as fin:
                generate(n,m,cap,cost,fin)
            p = subprocess.run(
                ['./benchmark/benchmark-mcf'],
                stdin=open(fname,'r'),capture_output=True)
            data = p.stdout.decode().split('\n')
            for d in data:
                try:
                    name,val = d.split()
                    val = int(val)
                    series[name][i] = max(series[name][i],val)
                except:
                    pass
    print(series)
    plot_data(series)
