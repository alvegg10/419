import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt


#   CPP SOLVER CLASS    ------------------------------------------------------------------------------------------------
class CPP:
    def __init__(self, vertices):
        self.N = vertices
        self.delta = []
        self.neg = []
        self.pos = []
        self.arcs = []
        self.label = []
        self.f = []  # repeated arcs
        self.c = []  # costs of cheapest paths
        self.cheapestLabel = []  # labels of cheapest arcs
        self.defined = []  # whether path cost is defined between vertices
        self.path = []  # spanning tree of the graph
        self.basicCost = 0.0  # total cost of traversing each arc once

        if self.N <= 0:
            raise ValueError("Graph is empty")

        self.delta = [0] * self.N
        self.defined = [[False] * self.N for _ in range(self.N)]
        self.label = [[[] for _ in range(self.N)] for _ in range(self.N)]
        self.c = [[0.0] * self.N for _ in range(self.N)]
        self.f = [[0] * self.N for _ in range(self.N)]
        self.arcs = [[0] * self.N for _ in range(self.N)]
        self.cheapestLabel = [[""] * self.N for _ in range(self.N)]
        self.path = [[0] * self.N for _ in range(self.N)]
        self.basicCost = 0
        self.NONE = -1  # anything < 0

    def addArc(self, lab, u, v, cost):
        if not self.defined[u][v]:
            self.label[u][v] = []
        self.label[u][v].append(lab)
        self.basicCost += cost

        if not self.defined[u][v] or self.c[u][v] > cost:
            self.c[u][v] = cost
            self.cheapestLabel[u][v] = lab
            self.defined[u][v] = True
            self.path[u][v] = v
            self.arcs[u][v] += 1
            self.delta[u] += 1
            self.delta[v] -= 1

        return self

    def solve(self):
        self.leastCostPaths()
        self.checkValid()
        self.findUnbalanced()
        self.findFeasible()
        while self.improvements():
            pass

    def leastCostPaths(self):
        for k in range(self.N):
            for i in range(self.N):
                if self.defined[i][k]:
                    for j in range(self.N):
                        if self.defined[k][j] and (
                                not self.defined[i][j] or self.c[i][j] > self.c[i][k] + self.c[k][j]):
                            self.path[i][j] = self.path[i][k]
                            self.c[i][j] = self.c[i][k] + self.c[k][j]
                            self.defined[i][j] = True

                            if i == j and self.c[i][j] < 0:
                                return

    def checkValid(self):
        for i in range(self.N):
            for j in range(self.N):
                if not self.defined[i][j]:
                    raise ValueError("Graph is not strongly connected")
                if i == j and self.c[i][j] < 0:
                    raise ValueError("Graph has a negative cycle")

    def findUnbalanced(self):
        nn, np = 0, 0
        for i in range(self.N):
            if self.delta[i] < 0:
                nn += 1
            elif self.delta[i] > 0:
                np += 1

        self.neg = [0] * nn
        self.pos = [0] * np
        nn = np = 0

        for i in range(self.N):
            if self.delta[i] < 0:
                self.neg[nn] = i
                nn += 1
            elif self.delta[i] > 0:
                self.pos[np] = i
                np += 1

    def findFeasible(self):
        delta = [0] * self.N
        for i in range(self.N):
            delta[i] = self.delta[i]

        for u in range(len(self.neg)):
            i = self.neg[u]
            for v in range(len(self.pos)):
                j = self.pos[v]
                self.f[i][j] = -delta[i] if -delta[i] < delta[j] else delta[j]
                delta[i] += self.f[i][j]
                delta[j] -= self.f[i][j]

    def improvements(self):
        residual = CPP(self.N)

        for u in range(len(self.neg)):
            i = self.neg[u]
            for v in range(len(self.pos)):
                j = self.pos[v]
                residual.addArc(None, i, j, self.c[i][j])
                if self.f[i][j] != 0:
                    residual.addArc(None, j, i, -self.c[i][j])

        residual.leastCostPaths()  # find a negative cycle

        for i in range(self.N):
            if residual.c[i][i] < 0:  # cancel the cycle (if any)
                k = 0
                u, v = 0, 0
                kunset = True
                u = i
                # find k to cancel
                while (u := residual.path[u][i]) != i:
                    v = residual.path[u][i]
                    if residual.c[u][v] < 0 and (kunset or k > self.f[v][u]):
                        k = self.f[v][u]
                        kunset = False

                u = i
                while (u := residual.path[u][i]) != i:
                    v = residual.path[u][i]
                    if residual.c[u][v] < 0:
                        self.f[v][u] -= k
                    else:
                        self.f[u][v] += k

                return True

        return False

    def cost(self):
        return self.basicCost + self.phi()

    def phi(self):
        phi = 0
        for i in range(self.N):
            for j in range(self.N):
                phi += self.c[i][j] * self.f[i][j]
        return phi

    def findPath(self, from_, f):
        for i in range(self.N):
            if f[from_][i] > 0:
                return i
        return self.NONE

    def printCPT(self, startVertex):
        v = startVertex

        # Delete next 7 lines to be faster, but non-reentrant
        arcs = [[0] * self.N for _ in range(self.N)]
        f = [[0] * self.N for _ in range(self.N)]
        for i in range(self.N):
            for j in range(self.N):
                arcs[i][j] = self.arcs[i][j]
                f[i][j] = self.f[i][j]

        while True:
            u = v
            v = self.findPath(u, f)
            if v != self.NONE:
                f[u][v] -= 1  # remove path
                # break down path into its arcs
                p = None
                while u != v:
                    p = self.path[u][v]
                    print("Take arc", self.cheapestLabel[u][p], "from", u, "to", p)
                    u = p
            else:
                bridgeVertex = self.path[u][startVertex]
                if arcs[u][bridgeVertex] == 0:
                    break  # finished if bridge already used
                v = bridgeVertex
                # find an unused arc, using bridge last
                for i in range(self.N):
                    if i != bridgeVertex and arcs[u][i] > 0:
                        v = i
                        break
                arcs[u][v] -= 1  # decrement count of parallel arcs
                print("Take arc", self.label[u][v][arcs[u][v]], "from", u, "to", v)  # use each arc label in turn


#   MAIN    ------------------------------------------------------------------------------------------------------------

edges = pd.read_csv('TBL_CZONE_10.csv', header=0)
edges.to_dict()
print(edges)

#   Assigning unique values to unique nodes
df = pd.DataFrame(edges)

combined_labels = set(df['fnode'].unique()) | set(df['endnode'].unique())

vertex_to_int = {label: idx for idx, label in enumerate(combined_labels)}
#   Return the max unique value + 1 as the number of verticies
num_of_verticies = max(vertex_to_int.values()) + 1
print("NUMBER OF VERTICIES = ", num_of_verticies)

rev_edges = {}
for column in ["fnode", "endnode"]:
    rev_edges[column] = {idx: vertex_to_int[value] for idx, value in df[column].items()}

for column in ["lngth", "depot"]:
    rev_edges[column] = {idx: value for idx, value in df[column].items()}

#   Implement CPP
M = CPP(num_of_verticies)

#   Add the arcs
for idx in rev_edges['fnode']:
    fnode_val = rev_edges['fnode'][idx]
    endnode_val = rev_edges['endnode'][idx]
    lngth_val = rev_edges['lngth'][idx]
    print("ADDING ARC: ", idx, fnode_val, endnode_val, lngth_val)
    M.addArc(idx, fnode_val, endnode_val, lngth_val)
    print("...SUCCESS")

#   Solve the CPP and output results
M.solve()
M.printCPT(4)
print('cost = ', str(M.cost()))