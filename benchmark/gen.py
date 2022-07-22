import networkx 
import random
import sys
import tempfile
import numpy
import subprocess
import matplotlib.pyplot as plt
import json

def generate(n,m,mcap,mwei,f):
    S=0
    T=1
    g=networkx.gnm_random_graph(n,m,directed=True)
    print(n,m,S,T,file=f)
    for a,b in g.edges():
        print(a,b,random.choice(range(mcap)),random.choice(range(mwei)),file=f)

def plot_data(n_array,series,outname):
    fig = plt.figure()
    ax = fig.subplots()
    for x in series:
        ax.loglog(n_array,series[x],'-o',label=x)
    ax.legend()
    ax.set_xlabel('N nodes')
    ax.set_ylabel('Time (micro seconds)')
    ax.grid()
    fig.savefig(outname)
    plt.show()

def save_data(n_array,series,outname):
    alldata={'N' : n_array}
    alldata['series']=dict()
    for name in series:
        alldata['series'][name]=[ int(x) for x in series[name] ]
    json.dump(alldata,open(outname,'w'))

def read_data(filename):
    alldata = json.load(open(filename,'r'))
    return alldata['N'],alldata['series']

def cat(fname):
    f = open(fname,'r')
    for l in f:
        print(l,end='')

def run_benchmark():
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
            # cat(fname)
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
    return n_array, series

if __name__ == "__main__":
    #n_array,series = read_data('latest.json')
    n_array,series = run_benchmark()
    save_data(n_array,series,'latest.json')
    plot_data(n_array,series,'latest.png')
