import networkx as nx

def compute_upper_bound(n, m, edges, k):

    # tạo đồ thị gốc
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)

    # tạo complete graph
    G2 = nx.complete_graph(n)
    E_c = list(nx.non_edges(G))
    edge_checked = set()
    clique_unchanged = 0
    clique_changed_flag = 0
    lst2 = []

    while E_c:
        q_old = nx.graph_clique_number(G2)

        e_incident = None
        for e in E_c:
            if e[0] not in edge_checked and e[1] not in edge_checked:
                e_incident = e
                edge_checked.add(e[0])
                edge_checked.add(e[1])
                break

        if e_incident is None:
            if len(edge_checked) != n:
                for e in E_c:
                    if e[0] in edge_checked and e[1] in edge_checked:
                        e_incident = e
                        break
            else:
                edge_checked = set()

        if e_incident is None:
            edge_checked = set()
            continue

        E_c.remove(e_incident)
        if G2.has_edge(*e_incident):
            G2.remove_edge(*e_incident)
        elif G2.has_edge(e_incident[1], e_incident[0]):
            G2.remove_edge(e_incident[1], e_incident[0])

        q_now = nx.graph_clique_number(G2)

        if q_now == q_old and q_now >= int(n / 2) and clique_changed_flag == 0:
            clique_unchanged = q_old
            clique_changed_flag = 1

        if q_now >= clique_unchanged:
            upper_bound = (k - 1) * q_now + n + 1 - k
        else:
            upper_bound = (k - 1) * clique_unchanged + n + 1 - k

        lst2.append(upper_bound)

    return lst2[-1] if lst2 else None

if __name__ == "__main__":
    n = int(input())
    m = int(input())
    edges = [tuple(map(int, input().split())) for _ in range(m)]
    k = int(input())
    ub = compute_upper_bound(n,m,edges,k)
    print(ub)