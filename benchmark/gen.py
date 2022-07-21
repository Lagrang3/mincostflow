import networkx 
import random
import sys

def generate(n,m,mcap,mwei):
    S=0
    T=1
    g=networkx.gnm_random_graph(n,m,directed=True)
    print(n,m,S,T)
    for a,b in g.edges():
        print(a,b,random.choice(range(mcap)),random.choice(range(mwei)))

if __name__ == "__main__":
    n,m,mcap,mcost = map(int,sys.argv[1:])
    generate(n,m,mcap,mcost)

# generate('linear-100.txt',100,300,20,20)
# generate('linear-1000.txt',1000,3000,200,200)
# generate('quad-100.txt',100,2000,20,20)
# generate('quad-200.txt',200,8000,200,200)
# generate('quad-400.txt',400,32000,200,200)
